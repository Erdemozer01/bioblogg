from django.views import generic
from blog.models import Posts
from dash import html, dcc
from django_plotly_dash import DjangoDash
import plotly.express as px
import pandas as pd
from dash import dash_table
from collections import OrderedDict
from django.contrib import messages
from django.shortcuts import redirect, get_object_or_404, reverse
from django.conf import settings
from accounts.models import Profile,SocialMedia, ContactModel
from django.http import HttpResponseRedirect
from blog.forms import ProfileModelForm, SocialMediaModelForm
from django.db.models import Q

class DashBoardView(generic.ListView):
    template_name = "dashboard/pages/index.html"
    model = Posts
    paginate_by = 10

    def get(self, request, *args, **kwargs):
        if not self.request.user.is_staff:
            messages.error(request, "Yetkili girişi yapınız!")
            return redirect('%s?next=/blog/' % (settings.LOGIN_URL))

        post_graph = DjangoDash('PostGraph')
        post_table = DjangoDash('PostTable')

        post_data = OrderedDict(
            [
                ("#", [post.id for post in Posts.objects.all()]),
                ("Tarih", [post.created.date() for post in Posts.objects.all()]),
                ("Başlık", [post.title for post in Posts.objects.all()]),
                ("Yazar", [post.author.username for post in Posts.objects.all()]),
                ("Kategori", [post.category.title for post in Posts.objects.all()]),
                ("Likes", [post.likes.count() for post in Posts.objects.all()]),
                ("Dislikes", [post.dislike.count() for post in Posts.objects.all()]),
            ]
        )

        post_table_df = pd.DataFrame(
            OrderedDict([(name, col_data) for (name, col_data) in post_data.items()])
        )

        post_table.layout = html.Div([

            html.Br(),

            dash_table.DataTable(
                data=post_table_df.to_dict('records'),
                columns=[{'id': c, 'name': c} for c in post_table_df.columns],
                style_table={'height': '300px', 'overflowY': 'auto'},
                style_cell={'textAlign': 'center'},
                filter_action="native",
                filter_options={"placeholder_text": "Ara", 'case': 'insensitive'},
                sort_action="native",
                sort_mode="multi",
                page_action="native",
                page_current=0,
                page_size=10,
                cell_selectable=True
            )
        ])

        fig = px.bar(
            post_data,
            x="Kategori",
            color="Kategori",
            barmode="group",
            title="Gönderi-Kategori",
            labels={'count': 'Gönderi Sayısı'},
            hover_data=['Tarih', 'Başlık', 'Likes'],
            text_auto=True,
        )

        post_graph.layout = html.Div(children=[
            dcc.Graph(
                id='example-graph',
                figure=fig
            )
        ])

        return super().get(request, *args, **kwargs)

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super().get_context_data(**kwargs)
        context['messages_list'] = ContactModel.objects.filter(author=self.request.user.pk)
        context['messages_count'] = ContactModel.objects.filter(author=self.request.user.pk).count()
        return context


class ProfilePostListView(generic.ListView):
    template_name = 'dashboard/pages/index.html'
    model = Posts

    def get(self, request, *args, **kwargs):
        if not self.request.user.is_staff:
            messages.error(request, "Yetkili girişi yapınız!")
            return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
        return super().get(request, *args, **kwargs)

    def get_queryset(self):
        return Posts.objects.filter(author=self.request.user.pk)

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super().get_context_data(**kwargs)
        context['messages_list'] = ContactModel.objects.filter(author=self.request.user.pk)
        context['messages_count'] = ContactModel.objects.filter(author=self.request.user.pk).count()
        return context


class ProfileListView(generic.ListView):
    template_name = 'dashboard/pages/index.html'
    model = Profile

    def get(self, request, *args, **kwargs):
        if not self.request.user.is_staff:
            messages.error(request, "Yetkili girişi yapınız!")
            return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
        return super().get(request, *args, **kwargs)

    def get_queryset(self):
        return Profile.objects.filter(user=self.request.user.pk)

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super().get_context_data(**kwargs)
        context['user_post_list'] = Posts.objects.filter(author=self.request.user)
        context['messages_list'] = ContactModel.objects.filter(author=self.request.user.pk)
        context['messages_count'] = ContactModel.objects.filter(author=self.request.user.pk).count()
        return context


class ProfileUpdateView(generic.UpdateView):
    template_name = 'dashboard/pages/update.html'
    model = Profile
    form_class = ProfileModelForm

    def get(self, request, *args, **kwargs):
        if not self.request.user.is_staff:
            messages.error(request, "Yetkili girişi yapınız!")
            return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
        return super().get(request, *args, **kwargs)

    def form_valid(self, form):
        instance = form.save(commit=False)
        instance.user_id = self.request.user.id
        instance.save()
        messages.success(self.request, "Profiliniz başarılı bir şekilde güncellendi!")
        return HttpResponseRedirect(self.request.META.get('HTTP_REFERER'))

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['messages_list'] = ContactModel.objects.filter(author=self.request.user.pk)
        context['messages_count'] = ContactModel.objects.filter(author=self.request.user.pk).count()
        return context


class MessagesListView(generic.ListView):
    template_name = 'dashboard/pages/index.html'
    model = ContactModel
    context_object_name = 'messages_list'

    def get_queryset(self):
        return ContactModel.objects.filter(author=self.request.user.pk)

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super().get_context_data(**kwargs)
        context['messages_count'] = ContactModel.objects.filter(author=self.request.user.pk).count()
        return context


class MessagesDetailView(generic.DetailView):
    template_name = 'dashboard/pages/messages.html'
    model = ContactModel

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super().get_context_data(**kwargs)
        context['messages_list'] = ContactModel.objects.filter(author=self.request.user.pk)
        context['messages_count'] = ContactModel.objects.filter(author=self.request.user.pk).count()
        return context


def mark_read_message(request, pk):
    message = get_object_or_404(ContactModel, pk=pk)
    message.is_read = True
    message.save()
    return HttpResponseRedirect(reverse('dashboard:messages_list', kwargs={
        'user': request.user.id
    }))


def mark_as_read_all(request):
    for message in ContactModel.objects.filter(author=request.user.pk).all():
        message.is_read = True
        message.save()
    return HttpResponseRedirect(reverse('dashboard:messages_list', kwargs={
        'user': request.user.id
    }))


def report_message(request, pk):
    message = get_object_or_404(ContactModel, pk=pk)
    message.is_report = True
    message.save()
    messages.success(request, "Mesajı rapor ettiniz!")
    return HttpResponseRedirect(reverse('dashboard:messages_list', kwargs={
        'user': request.user.id
    }))


def message_delete(request, pk):
    obj = get_object_or_404(ContactModel, pk=pk)
    obj.delete()
    messages.success(request, "Mesaj başarılı bir şekilde silindi!")
    return HttpResponseRedirect(reverse('dashboard:messages_list', kwargs={
        'user': request.user.id
    }))

class SocialMediaListView(generic.ListView):
    template_name = 'dashboard/pages/index.html'
    model = SocialMedia

    def get_queryset(self):
        return SocialMedia.objects.filter(user=self.request.user.pk)

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super().get_context_data(**kwargs)
        context['messages_list'] = ContactModel.objects.filter(author=self.request.user.pk)
        context['messages_count'] = ContactModel.objects.filter(author=self.request.user.pk).count()
        return context


class SocialMediaCreateView(generic.CreateView):
    template_name = 'dashboard/pages/create.html'
    model = SocialMedia
    form_class = SocialMediaModelForm

    def get(self, request, *args, **kwargs):
        if not self.request.user.is_staff:
            messages.error(request, "Yetkili girişi yapınız!")
            return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
        return super().get(request, *args, **kwargs)

    def get_queryset(self):
        return SocialMedia.objects.filter(user=self.request.user.id)

    def form_valid(self, form):
        instance = form.save(commit=False)
        instance.user_id = self.request.user.id
        if SocialMedia.objects.filter(user=self.request.user.id, social=instance.social).exists():
            messages.error(self.request, "Sosyal Medya bilgisi mevcut!")
            return HttpResponseRedirect(reverse('dashboard:social_dash', kwargs={
                'user': self.request.user
            }))
        instance.save()
        messages.success(self.request, "Sosyal Medya Bilgisi başarılı bir şekilde Eklendi!")
        return HttpResponseRedirect(reverse('dashboard:social_dash', kwargs={
            'user': self.request.user
        }))

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super().get_context_data(**kwargs)
        context['messages_list'] = ContactModel.objects.filter(author=self.request.user.pk)
        context['messages_count'] = ContactModel.objects.filter(author=self.request.user.pk).count()
        return context


class SocialMediaUpdateView(generic.UpdateView):
    template_name = 'dashboard/pages/update.html'
    model = SocialMedia
    form_class = SocialMediaModelForm

    def get(self, request, *args, **kwargs):
        if not self.request.user.is_staff:
            messages.error(request, "Yetkili girişi yapınız!")
            return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
        return super().get(request, *args, **kwargs)

    def get_queryset(self):
        return SocialMedia.objects.filter(user=self.request.user.id)

    def form_valid(self, form):
        instance = form.save(commit=False)
        instance.user_id = self.request.user.id
        instance.save()
        messages.success(self.request, "Sosyal Medya Bilgisi başarılı bir şekilde güncellendi!")
        return HttpResponseRedirect(reverse('dashboard:social_dash', kwargs={
            'user': self.request.user
        }))

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super().get_context_data(**kwargs)
        context['messages_list'] = ContactModel.objects.filter(author=self.request.user.pk)
        context['messages_count'] = ContactModel.objects.filter(author=self.request.user.pk).count()
        return context
