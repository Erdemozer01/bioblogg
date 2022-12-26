from django.views import generic
from accounts.models import Profile
from django.contrib.auth.models import User
from blog.models import Posts
from hitcount.models import Hit, HitCount
from django.contrib import messages
from django.shortcuts import redirect
from django.conf import settings
from blog.models import Comments
from django_plotly_dash import DjangoDash
from dash.dependencies import Input, Output
from dash import html, dcc, dash_table
import dash_bootstrap_components as dbc
import plotly.express as px
import pandas as pd
from collections import OrderedDict


class DashboardView(generic.ListView):
    template_name = "dashboard/pages/dashboard.html"
    model = Profile

    def get(self, request, *args, **kwargs):
        if not request.user.is_staff:
            messages.error(request, "Yetkili girişi yapınız!")
            return redirect('%s?next=/blog/' % (settings.LOGIN_URL))

        if self.request.user.is_superuser:
            user_agent = []
            users = []
            created = []

            for hit in Hit.objects.all():
                user_agent.append(hit.user_agent[111:])
                created.append(hit.created.date())
                users.append(str(hit.user))

            hits_table = DjangoDash('Hits')

            data = OrderedDict(
                [
                    ("Tarih", created),
                    ("Kullanıcılar", users),
                    ("İP", user_agent),
                ]
            )

            hits_data = pd.DataFrame(
                OrderedDict(
                    [(name, col_data) for (name, col_data) in data.items()]
                )
            )

            hits_table.layout = html.Div([
                dbc.Label('Kullanıcı - IP'),
                dash_table.DataTable(
                    data=hits_data.to_dict('records'),
                    columns=[{'id': c, 'name': c} for c in hits_data.columns],
                    style_table={'height': '250px', 'overflowY': 'auto'},
                    style_cell={'textAlign': 'center'},
                    filter_action="native",
                    filter_options={"placeholder_text": "Ara"},
                    sort_action="native",
                    sort_mode="multi",
                    page_action="native",
                    page_current=0,
                    page_size=10,
                )
            ])

        user_graph = DjangoDash("UserGraph")

        users_table = DjangoDash("UserTable")

        date_joined = []
        user_model = []
        first_name = []
        last_name = []

        for user in User.objects.all():
            user_model.append(user.username)
            first_name.append(user.first_name)
            last_name.append(user.last_name)
            date_joined.append(user.date_joined.date())

        users_data = OrderedDict(
            [
                ("Katılma Tarihi", date_joined),
                ("Kullanıcılar", user_model),
                ("Adı", first_name),
                ("Soyadı", last_name),
            ]
        )

        users_data_frame = pd.DataFrame(
            OrderedDict([(name, col_data) for (name, col_data) in users_data.items()])
        )

        users_table.layout = html.Div([
            dbc.Label(f'Kullanıcılar Tablosu, Toplam Kullanıcı Sayısı: {len(user_model)}'),
            dash_table.DataTable(
                data=users_data_frame.to_dict('records'),
                columns=[{'id': c, 'name': c} for c in users_data_frame.columns],
                style_table={'height': '250px', 'overflowY': 'auto'},
                style_cell={'textAlign': 'center'},
                filter_action="native",
                filter_options={"placeholder_text": "Ara"},
                sort_action="native",
                sort_mode="multi",
                page_action="native",
                page_current=0,
                page_size=10,
            )
        ], className="container")

        fig = px.line(data_frame=users_data_frame, x="Katılma Tarihi")

        user_graph.layout = html.Div([
            dbc.Label('Kullanıcılar'),
            dcc.Graph(figure=fig)
        ])

        post_table_app = DjangoDash('post_table')

        post_pk = []

        post_title = []

        post_author = []

        post_date = []

        post_category = []

        for post in Posts.objects.all().order_by('-id'):
            post_pk.append(post.pk)
            post_title.append(post.title)
            post_author.append(post.author.username)
            post_date.append(post.created.date())
            post_category.append(post.category.title)

        posts_data = OrderedDict(
            [
                ('Sayı', range(len(post_pk))),
                ("Tarih", post_date),
                ("Yazar", post_author),
                ("Gönderi", post_title),
                ("Kategori", post_category),
            ]
        )

        posts_table = pd.DataFrame(
            OrderedDict([(name, col_data) for (name, col_data) in posts_data.items()])
        )

        post_table_app.layout = html.Div([
            dbc.Label('Gönderiler'),
            dash_table.DataTable(
                id='datatable-interactivity',
                data=posts_table.to_dict('records'),
                columns=[{'id': c, 'name': c} for c in posts_table.columns],
                style_table={'height': '250px', 'overflowY': 'auto'},
                style_cell={'textAlign': 'center'},
                filter_action="native",
                filter_options={"placeholder_text": "Ara"},
                sort_action="native",
                sort_mode="multi",
                page_action="native",
                page_current=0,
                page_size=10,
            ), html.Div(id='datatable-interactivity-container')
        ], className="container")

        post_bar_fig = px.bar(data_frame=posts_table, x="Kategori", color='Kategori')

        post_graph_bar = DjangoDash('post_graph_bar')

        post_graph_bar.layout = html.Div([
            dbc.Label('Gönderiler'),
            dcc.Graph(figure=post_bar_fig)
        ])

        hits_table_app = DjangoDash('hits_table_son')

        post_hits_object_son = []
        post_hits_object_list = []

        hits_post_title = []
        hits_post_author = []
        hits_post_created = []
        hits_post_category = []

        for hits_count_obj in HitCount.objects.all():
            post_hits_object_son.append(hits_count_obj.hits)
            post_hits_object_list.append(hits_count_obj.object_pk)

        for post_hits in post_hits_object_list:
            for posts in Posts.objects.filter(id=post_hits):
                hits_post_title.append(posts.title)
                hits_post_author.append(posts.author)
                hits_post_created.append(posts.created)
                hits_post_category.append(posts.category)

        print(hits_post_title)

        hits_data_son = OrderedDict(
            [

                ("index:", range(len(post_hits_object_son))),
                ("Tarih:", hits_post_created),
                ("Başlık:", hits_post_title),
                ("Yazar:", hits_post_author),
                ("Okunma:", post_hits_object_son),

            ]
        )

        print(post_hits_object_son)

        hits_df_son = pd.DataFrame(OrderedDict([(name, col_data) for (name, col_data) in hits_data_son.items()]))

        hits_table_app.layout = html.Div([
            dbc.Label('Gönderi Okunma'),
            dash_table.DataTable(
                id='datatable-hits',
                data=posts_table.to_dict('records'),
                columns=[{'id': c, 'name': c} for c in hits_df_son.columns],
                style_table={'height': '250px', 'overflowY': 'auto'},
                style_cell={'textAlign': 'center'},
                filter_action="native",
                filter_options={"placeholder_text": "Ara"},
                sort_action="native",
                sort_mode="multi",
                page_action="native",
                page_current=0,
                page_size=10,
            )
        ])

        return super().get(request, *args, **kwargs)

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super().get_context_data(**kwargs)
        context['comments'] = Comments.objects.filter(post__author=self.request.user.pk, status="okunmadı").order_by(
            'commentator__posts__comment').order_by('-id')[:5]
        context['likes'] = Posts.objects.filter(author=self.request.user).order_by('author__post_like__likes').order_by(
            '-id')
        return context


class BlogDashBoardView(generic.ListView):
    template_name = "dashboard/pages/blog.html"
    model = Posts

    def get(self, request, *args, **kwargs):

        if not request.user.is_staff:
            messages.error(request, "Yetkili girişi yapınız!")
            return redirect('%s?next=/blog/' % (settings.LOGIN_URL))

        staff_user_post_title = []

        staff_user_post_date = []
        staff_user_category = []

        for staff in Posts.objects.filter(author=self.request.user).order_by('-id'):
            staff_user_post_title.append(staff.title)
            staff_user_category.append(staff.category.title)
            staff_user_post_date.append(staff.created.date())

        blog_app = DjangoDash(f"{self.request.user}_posts", suppress_callback_exceptions=True)

        posts_data = OrderedDict(
            [
                ("Yayınlanma Tarihi", staff_user_post_date),
                ("Gönderi Başlığı", staff_user_post_title),
                ("Kategori", staff_user_category),
            ]
        )

        posts_data_frame = pd.DataFrame(posts_data)

        blog_app.layout = html.Div([
            dbc.Label('Gönderilerim'),
            dash_table.DataTable(
                data=posts_data_frame.to_dict('records'),
                columns=[{'id': c, 'name': c} for c in posts_data_frame.columns],
                style_table={'height': '250px', 'overflowY': 'auto'},
                style_cell={'textAlign': 'center'},
                filter_action="native",
                filter_options={"placeholder_text": "Ara"},
                sort_action="native",
                sort_mode="multi",
                page_action="native",
                page_current=0,
                page_size=10,
            )
        ])

        return super().get(request, *args, **kwargs)

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super().get_context_data(**kwargs)
        context['comments'] = Comments.objects.filter(post__author=self.request.user.pk, status="okunmadı").order_by(
            'commentator__posts__comment').order_by('-id')[:5]
        context['likes'] = Posts.objects.filter(author=self.request.user).order_by('author__post_like__likes').order_by(
            '-id')
        context['staff_posts'] = Posts.objects.filter(author=self.request.user.is_staff)
        context['posts'] = f"{self.request.user}_posts"
        return context
